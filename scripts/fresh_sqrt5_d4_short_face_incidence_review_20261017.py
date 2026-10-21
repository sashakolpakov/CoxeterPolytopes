#!/usr/bin/env python3
"""Independent exact reconstruction of the D4 incidence argument's local inputs."""
import argparse
import json
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_d4_short_face_incidence_review_20261017'
BASE = 'fresh_sqrt5_d4_short_face_incidence_20261017'
PHI = (1 + s.sqrt(5)) / 2


def inertia(matrix):
    """Exact congruence elimination, including zero diagonal pivots."""
    m = s.Matrix(matrix)
    counts = [0, 0, 0]
    while m.rows:
        pivot = next((i for i in range(m.rows) if m[i, i] != 0), None)
        if pivot is None:
            pair = next(((i, j) for i, j in combinations(range(m.rows), 2)
                         if m[i, j] != 0), None)
            if pair is None:
                counts[2] += m.rows
                break
            i, j = pair
            change = s.eye(m.rows)
            change[j, i] = 1
            m = change.T * m * change
            pivot = i
        ids = [pivot] + [i for i in range(m.rows) if i != pivot]
        m = m.extract(ids, ids)
        sign = s.sign(s.simplify(m[0, 0]))
        assert sign in (-1, 1)
        counts[0 if sign == 1 else 1] += 1
        col = m[1:, 0]
        m = (m[1:, 1:] - col * col.T / m[0, 0]).applyfunc(s.simplify)
    return counts


def d4():
    g = 2 * s.eye(4)
    for j in (1, 2, 3):
        g[0, j] = g[j, 0] = -1
    return g


def extend(g, row):
    return g.row_join(s.Matrix(row)).col_join(s.Matrix([list(row) + [2]]))


def run():
    g = d4()
    assert inertia(g) == [4, 0, 0]
    attachments = []
    for port in range(4):
        keep = [i for i in range(4) if i != port]
        matches = []
        for row in product((0, -1), repeat=3):
            e = extend(g.extract(keep, keep), row)
            degrees = sorted(sum(e[i, j] != 0 for j in range(4) if j != i)
                             for i in range(4))
            if degrees == [1, 1, 1, 3]:
                matches.append(list(row))
        assert matches == [[int(g[port, i]) for i in keep]]
        attachments.append(dict(port=port, unique_contacts=matches[0]))

    triangles = []
    survivors = []
    for port in range(4):
        for label, value in ((2, s.Integer(0)), (3, s.Integer(-1)), (5, -PHI)):
            row = [value if j == port else g[port, j] for j in range(4)]
            big = extend(g, row)
            physical = inertia(big)
            conjugate = inertia(big.xreplace({s.sqrt(5): -s.sqrt(5)}))
            for other in range(4):
                if other == port:
                    continue
                keep = [i for i in range(5) if i != other]
                third = inertia(big.extract(keep, keep))
                accepted = physical == [4, 1, 0] and conjugate == [5, 0, 0] and third == [4, 0, 0]
                record = dict(port=port, other=other, label=label,
                              physical=physical, conjugate=conjugate,
                              third_corner=third, accepted=accepted)
                triangles.append(record)
                if accepted:
                    survivors.append([port, other, label])
    assert survivors == [[1, 0, 5], [2, 0, 5], [3, 0, 5]]
    assert inertia(s.Matrix([[2, -1, -PHI], [-1, 2, -1], [-PHI, -1, 2]])) == [2, 1, 0]

    a, b, c, li, lj = s.symbols('a b c lambda_i lambda_j')
    residual = s.Matrix([[a, c, a-li, c], [c, b, c, b-lj],
                         [a-li, c, a, c], [c, b-lj, c, b]])
    # Pair differences first, then pair sums. This displays the signature
    # information that is not supplied by the determinant identity alone.
    change = s.Matrix([[1, 0, 1, 0], [0, 1, 0, 1],
                       [-1, 0, 1, 0], [0, -1, 0, 1]])
    congruent = (change.T * residual * change / 2).applyfunc(s.expand)
    expected = s.diag(li, lj, 0, 0)
    expected[2:, 2:] = s.Matrix([[2*a-li, 2*c], [2*c, 2*b-lj]])
    assert congruent == expected
    determinant = li*lj*((2*a-li)*(2*b-lj)-4*c*c)
    assert s.expand(residual.det()-determinant) == 0
    schur = []
    for copies in (1, 2):
        big = s.diag(*([g] * copies))
        inv = big.inv()
        pairs = combinations(range(4), 2) if copies == 1 else product(range(4), range(4, 8))
        for i, j in pairs:
            h = inv.extract([i, j], [i, j]).inv()
            if copies == 2:
                assert h[0, 1] == 0
            elif i != 0:
                assert h == s.Matrix([[s.Rational(4, 3), -s.Rational(2, 3)],
                                      [-s.Rational(2, 3), s.Rational(4, 3)]])
            schur.append(dict(copies=copies, ports=[i, j], matrix=[[str(x) for x in h.row(k)] for k in range(2)]))
    assert len(schur) == 22
    gi, gj = s.symbols('g_i g_j')
    leaf_equation = ((2*a-li)*(2*b-lj)-4*c*c).subs(
        {a:s.Rational(4,3), b:s.Rational(4,3), c:-s.Rational(2,3), li:2-gi, lj:2-gj})
    assert s.expand(9*leaf_equation - ((-3*gi-2)*(-3*gj-2)-16)) == 0

    norm = lambda x, y: x*x+x*y-y*y
    assert [(x,y) for x,y in product(range(2),repeat=2) if norm(x,y)%2 == 0] == [(0,0)]
    x, y = s.symbols('x y', integer=True)
    residue = s.Poly(s.expand(norm(1+3*x, 3*y)-1), x, y)
    assert all(int(v)%3 == 0 for v in residue.coeffs())
    negative_norms = [-4**j for j in range(5)]
    assert all(256 % abs(n) == 0 and n % 3 == 2 for n in negative_norms)

    k = s.symbols('k', integer=True, positive=True)
    corners = 4*k*(4*k-1)/2
    threshold = s.Rational(5,2)*(4*k*k-2*k)-corners
    special_cost = 3*k/s.Integer(2)
    assert s.expand(threshold) == 2*k*k-3*k
    assert s.simplify(2*(threshold-special_cost)/k) == 4*k-9
    assert s.simplify(2*(corners-threshold)/k) == 12*k+2
    paths = [f'scripts/{TAG}.py', f'results/{TAG}.md']
    paths += [f'{folder}/{BASE}.{ext}' for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    paths += ['results/fresh_sqrt5_discount_structural_floor_20261016.md']
    return dict(status='PASS_INDEPENDENT_D4_RECONSTRUCTION_AND_WRITTEN_INCIDENCE_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                reattachments=attachments, triangle_cases=triangles,
                allowed_triangle_ports=survivors, schur_cases=schur,
                residual_congruence=[[str(t) for t in congruent.row(i)] for i in range(4)],
                negative_norm_options=negative_norms,
                frequency_inequality='(4*k-9)*V_D < (12*k+2)*V_other for k>=3',
                scope='Independent local arithmetic and written geometric/counting review. Published arithmeticity, finite Coxeter classification and simple-polytope h-vector facts are mathematical inputs. Not a formal verification or general dimension theorem.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], flush=True)
