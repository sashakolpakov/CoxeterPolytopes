#!/usr/bin/env python3
"""Independent E7/E8 finite-contact census over Q(sqrt2).

No producer imports. SymPy algebraic-field dot products reconstruct costs;
rational interval signs replace the producer's squared-sign comparison.
Every accepted replacement is also checked by a direct Gram determinant.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_arithmetic_sqrt2_local_20261015_review'
PRIMARY = 'fresh_arithmetic_sqrt2_local_20261015'
K = s.QQ.algebraic_field(s.sqrt(2))


def rational(x):
    return F(int(x.numerator), int(x.denominator))


def coefficients(x):
    row = [rational(t) for t in x.to_list()]
    row = [F(0)]*(2-len(row))+row
    return row[1], row[0]


def encoded(x):
    return list(map(str, coefficients(x)))


SIGN_CACHE = {}
def exact_sign(x, conjugate=False):
    a, b = coefficients(x)
    if conjugate:
        b = -b
    key = (a, b)
    if key in SIGN_CACHE:
        return SIGN_CACHE[key]
    if not b:
        answer = (a > 0)-(a < 0)
    else:
        lo, hi = F(7, 5), F(3, 2)
        assert lo*lo < 2 < hi*hi
        for _ in range(256):
            left, right = sorted((a+b*lo, a+b*hi))
            if left > 0:
                answer = 1
                break
            if right < 0:
                answer = -1
                break
            mid = (lo+hi)/2
            if mid*mid < 2:
                lo = mid
            else:
                hi = mid
        else:
            raise ArithmeticError(key)
    SIGN_CACHE[key] = answer
    return answer


def matrix_for_rank(n):
    # Independent arm construction followed by the literal producer ordering:
    # centre2, left arm1,0; right arm3,...,n-2; short armn-1.
    arms = [[2, 1, 0], [2]+list(range(3, n-1)), [2, n-1]]
    G = 2*s.eye(n)
    for arm in arms:
        for a, b in zip(arm, arm[1:]):
            G[a, b] = G[b, a] = -1
    assert G.det() == (2 if n == 7 else 1)
    assert all(G[:j, :j].det() > 0 for j in range(1, n+1))
    return G


def field_quad(inverse, contacts):
    value = K.zero
    active = [i for i, x in enumerate(contacts) if x]
    for i in active:
        for j in active:
            value += contacts[i]*inverse[i][j]*contacts[j]
    return value


def census(n, primary):
    G = matrix_for_rank(n)
    J = G.inv(method='DM')
    assert G.tolist() == primary['gram']
    assert [[str(J[i, j]) for j in range(n)] for i in range(n)] == primary['inverse']
    inverse = [[K.convert(J[i, j]) for j in range(n)] for i in range(n)]
    counts, accepted, direct_records, subcritical = {}, [], [], []
    covered = {p: [] for p in combinations(range(n), 2)}
    for kind, alphabet in [('mixed34', (0, 1, 2)), ('pure6', (0, 1)), ('pure8', (0, 1))]:
        counts[kind] = dict(total=0, hyperbolic=0, conjugate_positive=0,
                            with_two_finite_replacements=0)
        for row in product(alphabet, repeat=n):
            counts[kind]['total'] += 1
            if kind == 'mixed34':
                labels = [2 if z == 0 else 3 if z == 1 else 4 for z in row]
                contacts = [K.zero if z == 0 else K.one if z == 1 else K.unit for z in row]
                cost = field_quad(inverse, contacts)
                physical = [s.Integer(0) if z == 0 else s.Integer(1) if z == 1 else s.sqrt(2) for z in row]
            else:
                labels = [int(kind[4:]) if z else 2 for z in row]
                binary = [K.one if z else K.zero for z in row]
                unit = K.convert(3 if kind == 'pure6' else 2+s.sqrt(2))
                cost = unit*field_quad(inverse, binary)
                physical = [s.sqrt(3) if kind == 'pure6' else s.sqrt(2+s.sqrt(2)) for _ in row]
                physical = [z*x for z, x in zip(row, physical)]
            hyp = exact_sign(cost-K.convert(2))
            if hyp <= 0:
                if any(row):
                    subcritical.append(dict(kind=kind, labels=labels, cost=encoded(cost),
                                            relation='equal' if hyp == 0 else 'below'))
                continue
            counts[kind]['hyperbolic'] += 1
            if exact_sign(K.convert(2)-cost, conjugate=True) <= 0:
                continue
            counts[kind]['conjugate_positive'] += 1
            # Direct full Gram and direct replacement determinants: no
            # inverse-column replacement formula from the producer is used.
            extension = G.row_join(-s.Matrix(physical))
            extension = extension.col_join(s.Matrix([[-x for x in physical]+[2]]))
            full_det = s.expand(extension.det(method='domain-ge'))
            assert K.convert(full_det) == K.convert(G.det())*(K.convert(2)-cost)
            finite, replacement, det_rows = [], {}, []
            for i in range(n):
                old = [j for j in range(n) if j != i]
                selected = old+[n]
                determinant = s.expand(extension.extract(selected, selected).det(method='domain-ge'))
                old_determinant = G.extract(old, old).det()
                schur = K.convert(s.expand(determinant/old_determinant))
                if exact_sign(schur) > 0:
                    finite.append(i)
                replacement[str(i)] = encoded(schur)
                det_rows.append(dict(removed=i, determinant=encoded(K.convert(determinant)),
                                     old_determinant=str(old_determinant)))
            rec = dict(kind=kind, labels=labels, projection_cost=encoded(cost),
                       finite_replacements=finite, replacement_schur=replacement)
            index = len(accepted)
            accepted.append(rec)
            direct_records.append(dict(profile=index, full_determinant=encoded(K.convert(full_det)),
                                       replacements=det_rows))
            if len(finite) >= 2:
                counts[kind]['with_two_finite_replacements'] += 1
            for pair in combinations(finite, 2):
                covered[pair].append(index)
    assert counts == primary['counts']
    assert accepted == primary['surviving_profiles']
    absent = [list(p) for p, values in covered.items() if not values]
    assert absent == primary['pairs_without_triangle']
    assert sum(bool(v) for v in covered.values()) == primary['pairs_with_triangle']
    for r in primary['marked_pairs']:
        assert covered[tuple(r['pair'])] == r['profiles']
    # Complete scope check for the orthogonal-component loophole.
    expected = [(5, 'below', ['3/2', '0']), (0, 'equal', ['2', '0'])] if n == 7 else [(6, 'equal', ['2', '0'])]
    actual = []
    for row in subcritical:
        assert row['kind'] == 'mixed34' and row['labels'].count(3) == 1 and row['labels'].count(2) == n-1
        actual.append((row['labels'].index(3), row['relation'], row['cost']))
    assert sorted(actual) == sorted(expected)
    return dict(type=f'E{n}', counts=counts, profile_count=len(accepted),
                covered_pairs=[list(p) for p, values in covered.items() if values],
                excluded_pairs=absent, subcritical_nonzero_profiles=subcritical,
                direct_gram_determinants=direct_records,
                inverse_diagonal=[str(J[i, i]) for i in range(n)],
                inverse_strictly_positive=all(z > 0 for z in J))


def field_guards():
    # phi(m)^2/m is the product of (p-1)^2*p^(e-2). Only p=2,e=1
    # has a factor below1, namely1/2. Thus phi(m)<=4 implies m<=32.
    assert F((2-1)**2, 2) == F(1, 2)
    assert F((3-1)**2, 3) > 1
    # For p>=3, (p-1)^2-p=p(p-3)+1>0; increasing e multiplies by p.
    candidates = [m for m in range(3, 33) if sum(gcd(i, m) == 1 for i in range(1, m)) <= 4]
    assert candidates == [3, 4, 5, 6, 8, 10, 12]
    square_norms = [2, 8, 18]
    assert all(isqrt(z)**2 != z for z in square_norms)
    for value in (F(3), F(3, 2), F(6)):
        assert isqrt(value.numerator)**2 != value.numerator or isqrt(value.denominator)**2 != value.denominator
    # A high6 bond adjoining a nonright3 already gives a singular path;
    # label8 increases the first squared cosine and makes it indefinite.
    path6 = s.Matrix([[1, -s.sqrt(3)/2, 0], [-s.sqrt(3)/2, 1, -s.Rational(1, 2)], [0, -s.Rational(1, 2), 1]])
    assert path6.det() == 0
    path8 = path6.subs(s.sqrt(3), s.sqrt(2+s.sqrt(2)))
    assert s.simplify(path8.det()-(1-s.sqrt(2))/4) == 0
    return dict(phi_candidates=candidates, phi_bound='phi(m)^2>=m/2, hence m<=32',
                crossclass_square_norms=square_norms,
                high_path_determinants=['0', '(1-sqrt(2))/4'],
                ambient_extension_review='E7 only subcritical singleton3 at node5; extra attached finite trees have rational Schur costs. E8 has no nonzero subcritical profile.')


def run():
    data = json.loads((ROOT/'results'/f'{PRIMARY}.json').read_text())
    records = [census(n, p) for n, p in zip((7, 8), data['records'])]
    paths = [f'scripts/{PRIMARY}.py', f'results/{PRIMARY}.json', f'results/{PRIMARY}.md', f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_SQRT2_E7_E8_TRIANGLE_AND_AMBIENT_SCOPE_REVIEW',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                arithmetic='SymPy QQ(sqrt2) dot products and direct Gram determinants; exact rational interval signs; no producer imports.',
                field_guards=field_guards(), records=records,
                totals=dict(profiles=sum(sum(x['total'] for x in r['counts'].values()) for r in records),
                            accepted=sum(r['profile_count'] for r in records),
                            direct_replacement_grams=sum(len(x['replacements']) for r in records for x in r['direct_gram_determinants'])),
                scope='Complete necessary finite-contact profiles, direct local-triangle checks, and proof review of the orthogonal-component extension; no quadrilateral or dimension theorem.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(out.read_text()) == result, 'saved review mismatch'
        print('PASS saved independent replay')
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
        print(result['status'])
    print(json.dumps(result['totals']))
