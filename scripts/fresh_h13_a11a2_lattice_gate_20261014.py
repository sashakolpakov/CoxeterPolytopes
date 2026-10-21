#!/usr/bin/env python3
"""All integral overlattices of the normalized A11+A2 face lattice.

Exact finite discriminant enumeration, with an independent ambient-coordinate
realization of the only nontrivial integral extension. This is a local lattice
gate; it says nothing about reflectivity or arbitrary geometric completions.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h13_a11a2_lattice_gate_20261014'
SOURCE = 'results/fresh_h20_e7_face_20261014.json'


def dot(a, b):
    return sum(x*y for x, y in zip(a[:-1], b[:-1])) - 2*a[-1]*b[-1]


def determinant(a):
    a = [list(map(F, r)) for r in a]
    d = F(1)
    for i in range(len(a)):
        j = next(j for j in range(i, len(a)) if a[j][i])
        if i != j:
            a[i], a[j] = a[j], a[i]
            d = -d
        pivot = a[i][i]
        d *= pivot
        for j in range(i+1, len(a)):
            f = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= f*a[i][k]
    return d


def run():
    source = json.loads((ROOT/SOURCE).read_text())
    Q = [[F(x) for x in r] for r in source['intrinsic_normalization']['form']]
    # Independent Euclidean sum-zero models for A11 and A2, plus norm -2.
    basis = []
    for start, n in [(0, 12), (12, 3)]:
        for i in range(start, start+n-1):
            v = [F(0)]*16
            v[i], v[i+1] = F(1), F(-1)
            basis.append(v)
    basis.append([F(0)]*15+[F(1)])
    assert [[dot(a, b) for b in basis] for a in basis] == Q
    assert determinant(Q) == -72
    # These dual generators have exact orders 12, 3, 2 modulo the root lattice.
    g = [[F(11, 12)]+[F(-1, 12)]*11+[F(0)]*4,
         [F(0)]*12+[F(2, 3), F(-1, 3), F(-1, 3), F(0)],
         [F(0)]*15+[F(1, 2)]]
    assert all(dot(v, b).denominator == 1 for v in g for b in basis)
    assert [[dot(a, b) for b in g] for a in g] == [
        [F(11, 12), F(0), F(0)], [F(0), F(2, 3), F(0)],
        [F(0), F(0), F(-1, 2)]]
    classes = []
    integral = []
    even = []
    for a, b, c in product(range(12), range(3), range(2)):
        v = [a*g[0][i]+b*g[1][i]+c*g[2][i] for i in range(16)]
        q = dot(v, v)
        assert q == F(11*a*a, 12)+F(2*b*b, 3)-F(c*c, 2)
        classes.append(dict(coset=[a, b, c], norm_mod_two=str(q % 2)))
        if q.denominator == 1:
            integral.append([a, b, c])
            if q % 2 == 0:
                even.append([a, b, c])
    assert integral == [[0, 0, 0], [6, 0, 0]]
    assert even == [[0, 0, 0]]
    # omega_6 is in class 6*g[0]; its norm is 3, not 33.
    omega6 = [F(1, 2)]*6+[F(-1, 2)]*6+[F(0)]*4
    delta = [omega6[i]-6*g[0][i] for i in range(16)]
    assert all(x.denominator == 1 for x in delta)
    assert sum(delta[:12]) == sum(delta[12:15]) == delta[15] == 0
    assert dot(omega6, omega6) == 3
    extended = [omega6]+basis[1:]
    gram = [[dot(a, b) for b in extended] for a in extended]
    assert all(x.denominator == 1 for row in gram for x in row)
    assert determinant(gram) == -18
    # The norm-six mirror alpha-beta swaps the two A2 simple roots. It
    # preserves both lattices and splits the original A2 chamber in half.
    alpha, beta = basis[11:13]
    long = [a-b for a, b in zip(alpha, beta)]
    assert dot(long, long) == 6
    def reflect(v):
        return [x-dot(v, long)*r/3 for x, r in zip(v, long)]
    assert reflect(alpha) == beta and reflect(beta) == alpha
    assert all(reflect(v) == v for i, v in enumerate(basis) if i not in [11, 12])
    assert reflect(omega6) == omega6
    left = [F(0)]*12+[F(-4, 3), F(-1, 3), F(5, 3), F(0)]
    right = [F(0)]*12+[F(-5, 3), F(1, 3), F(4, 3), F(0)]
    assert all(dot(v, r) < 0 for v in [left, right] for r in [alpha, beta])
    assert dot(left, long) == 1 and dot(right, long) == -1
    paths = [SOURCE, f'scripts/{TAG}.py']
    return dict(status='PASS_EXACT_INTEGRAL_OVERLATTICE_AND_APEX_MIRROR_GATE',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        lattice='A11 + A2 + <-2>', determinant=-72,
        discriminant_orders=[12, 3, 2], discriminant_classes=classes,
        integral_self_pairing_classes=integral, even_classes=even,
        integral_overlattices=[dict(index=1, determinant=-72, even=True),
            dict(index=2, determinant=-18, even=False,
                 generator_ambient=list(map(str, omega6)),
                 basis_gram=[[str(x) for x in row] for row in gram])],
        norm_six_mirror=list(map(str, long)),
        mirror_action='swaps the two A2 roots; fixes A11 and time; preserves both integral overlattices',
        conclusion='No integral overlattice of this fixed normalized lattice removes the additional A2-bisecting reflection. Full reflection groups of these lattices do not retain the A11+A2 apex chamber.',
        scope='Does not exclude a restricted reflection subgroup, a different or rescaled lattice, algebraic source transfer, or any arbitrary geometric completion.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
