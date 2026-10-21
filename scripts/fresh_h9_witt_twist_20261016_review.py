#!/usr/bin/env python3
"""Independent local invariant, norm pool and zero-sector checks for new H9 host."""
import argparse
import json
import sys
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from coxeter_polytopes.exact.golden import Golden as K

TAG = 'fresh_h9_witt_twist_20261016_review'
PRIMARY = 'results/fresh_h9_witt_twist_20261016.json'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def enc(x):
    if isinstance(x, K):
        return x.to_json()
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [enc(v) for v in x]
    return x


def dec(v):
    return [K.from_json(x) for x in v]


def vector_key(v):
    return tuple((x.a, x.b) for x in v)


def inverse(matrix):
    n = len(matrix)
    a = [[Q(x) for x in row]+[Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    d = Q(1)
    for j in range(n):
        i = next(i for i in range(j, n) if a[i][j])
        if i != j:
            a[i], a[j] = a[j], a[i]
            d = -d
        p = a[j][j]
        d *= p
        a[j] = [x/p for x in a[j]]
        for i in range(n):
            if i != j:
                p = a[i][j]
                a[i] = [x-p*y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a], d


def zero_block(n):
    # q=sum(a_i²+b_i²)<=2: enumerate all sparse coefficient vectors,
    # independently of the primary dynamic count and known Bn root list.
    roots, count = [], 0
    for q in (1, 2):
        for support in combinations(range(2*n), q):
            for signs in product((-1, 1), repeat=q):
                coefficients = [0]*(2*n)
                for i, s in zip(support, signs):
                    coefficients[i] = s
                count += 1
                v = [K(coefficients[i], coefficients[n+i]) for i in range(n)]
                if sum((x*x for x in v), K()) == q:
                    assert all(x.b == 0 for x in v)
                    roots.append(v)
    assert len(roots) == 2*n*n
    return roots, count


def run():
    primary = json.loads((ROOT / PRIMARY).read_text())
    assert primary['status'] == 'PASS_DISTINCT_LOCAL_WITT_H9_HOST_COMPLETE_NORMS_AND_ZERO_SECTOR'
    for path, digest in primary['sources'].items():
        assert sha(path) == digest, path
    phi = K(0, 1)
    beta = 4-phi
    assert phi*phi == phi+1 and phi.norm() == -1
    assert beta.norm() == 11 and beta > 0 and beta.conjugate() > 0
    diagonal = [-phi]+[K(1)]*7+[beta]*2
    form = [[diagonal[i] if i == j else K() for j in range(10)] for i in range(10)]
    assert form == [dec(row) for row in primary['form']]
    assert diagonal[0] < 0 and all(x > 0 for x in diagonal[1:])
    assert all(x.conjugate() > 0 for x in diagonal)
    determinant = K(1)
    for d in diagonal:
        determinant *= d
    assert determinant == -phi*beta*beta and determinant.norm() == -121
    # beta*G^-1 integral implies q divides2beta for every primitive reflector.
    dual_multiplier = [beta/d for d in diagonal]
    assert all(x.is_integral for x in dual_multiplier)
    norms = [K(1), K(2), beta, 2*beta]
    assert norms == dec(primary['norm_and_zero']['primitive_norm_classes'])
    assert all(x > 0 and x.conjugate() > 0 for x in norms)
    assert all((2*beta/x).is_integral for x in norms)
    assert all((r*r-r-1) % 2 for r in range(2))  # 2 inert, hence prime.

    # Residue map O/(beta)=F11 sends phi to4. Two independent Hensel lifts:
    # phi is a local square; I4 is isotropic and, by determinant, split.
    assert (4*4-4-1) % 11 == 0 and (2*4-1) % 11 != 0
    assert (2*2-4) % 11 == 0 and 2*2 % 11 != 0
    binary_zeros = [(a, b) for a, b in product(range(11), repeat=2)
                    if (a*a+b*b) % 11 == 0]
    assert binary_zeros == [(0, 0)]
    seed = (1, 1, 3, 0)
    assert sum(a*a for a in seed) % 11 == 0 and 2*seed[0] % 11 != 0
    assert -1 % 11 not in {a*a % 11 for a in range(11)}
    assert Q(4, 11) < 1

    # Explicit E8-to-I8 rational isometry supplied as a witness, fully checked.
    basis = [[Q(x) for x in row] for row in primary['family_comparison']['E8_euclidean_simple_roots']]
    _, basis_det = inverse(basis)
    gram = [[sum(x*y for x, y in zip(a, b)) for b in basis] for a in basis]
    edges = {(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)}
    expected = [[Q(2) if i == j else Q(-1) if (min(i, j), max(i, j)) in edges else Q()
                 for j in range(8)] for i in range(8)]
    assert gram == expected and abs(basis_det) == 1
    inv, d = inverse(gram)
    assert d == 1 and all(x.denominator == 1 for row in inv for x in row)

    left, left_count = zero_block(7)
    right, right_count = zero_block(2)
    roots = [[K()]+r+[K()]*2 for r in left]+[[K()]*8+r for r in right]
    assert len(roots) == 106
    assert {vector_key(r) for r in roots} == {vector_key(dec(r)) for r in primary['norm_and_zero']['zero_roots']}

    def pair(a, b):
        return sum((d*x*y for d, x, y in zip(diagonal, a, b)), K())

    counts = Counter()
    for r in roots:
        q = pair(r, r)
        counts[norms.index(q)] += 1
        assert all((2*d*x/q).is_integral for d, x in zip(diagonal, r))
    assert [counts[i] for i in range(4)] == [14, 84, 4, 4]
    simple = []
    for start, n in ((1, 7), (8, 2)):
        for j in range(n-1):
            r = [K()]*10
            r[start+j], r[start+j+1] = K(-1), K(1)
            simple.append(r)
        r = [K()]*10
        r[start+n-1] = K(-1)
        simple.append(r)
    assert simple == [dec(r) for r in primary['norm_and_zero']['zero_simple_walls']]
    simple_gram = [[pair(a, b) for b in simple] for a in simple]
    labels = []
    for i, j in combinations(range(9), 2):
        assert simple_gram[i][j] <= 0
        squared = 4*simple_gram[i][j]*simple_gram[i][j]/(simple_gram[i][i]*simple_gram[j][j])
        labels.append(dict(pair=[i, j], label=2 if squared == 0 else 3 if squared == 1 else 4 if squared == 2 else None))
    assert all(r['label'] is not None for r in labels)
    assert labels == primary['norm_and_zero']['zero_relations']
    paths = [f'scripts/{TAG}.py', 'src/coxeter_polytopes/exact/golden.py', PRIMARY,
             'results/fresh_h9_witt_twist_20261016.md']
    return enc(dict(status='PASS_INDEPENDENT_DISTINCT_WITT_HOST_NORMS_AND_ZERO_SECTOR',
                    sources={p: sha(p) for p in paths}, form=form, determinant=determinant,
                    beta=beta, integral_dual_multiplier=dual_multiplier, primitive_norm_classes=norms,
                    residue_binary_zeros=binary_zeros, Hensel_I4_seed=seed,
                    E8_euclidean_basis=basis, E8_basis_determinant=basis_det,
                    local_Witt_index=3, E8_binary_local_Witt_index_at_least=4,
                    sparse_coefficient_cases=left_count+right_count,
                    zero_roots=roots, zero_counts=[counts[i] for i in range(4)],
                    zero_simple_walls=simple, zero_relations=labels,
                    scope='Independent exact arithmetic and finite residue/coefficient checks with the accompanying all-height norm and local-field proof. No positive-height enumeration or reflectivity/finite-volume conclusion. Local index separates this rational quadratic space from every E8-plus-binary form, including scalar similarities.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['local_Witt_index'], result['zero_counts'])
