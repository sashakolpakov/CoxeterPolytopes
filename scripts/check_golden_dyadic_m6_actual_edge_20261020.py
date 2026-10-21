#!/usr/bin/env python3
"""One exact finite M6 edge; run only through the parent-lock bounded wrapper.

No height prefix, global chamber enumeration, or compact H9 claim.
--write refuses replacement; --verify recomputes in a separate process.
"""
import argparse
import hashlib
import json
import signal
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from check_arithmetic_dn_opposed_fork_models_20261019 import field, inverse

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/golden_dyadic_m6_actual_edge_20261020.json'
PINS = {
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
        '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'results/golden_dyadic_m6_actual_short_edge_preparation_20261020.md':
        '2afd545c3c3166cd07381949b972fb012280b18ed673a090025b6aac30537dd2',
    'results/golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md':
        'd6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501',
    'results/p6_10_bugaenko_double_lineage_and_classification_audit.md':
        '651176b56e7d0d0a5e7596b1c889e187f0f5a9703b8336c456aa2e0b445d56c3',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def rational_sqrt(x):
    if x < 0:
        return None
    n, d = isqrt(x.numerator), isqrt(x.denominator)
    return F(n, d) if n*n == x.numerator and d*d == x.denominator else None


def square_roots(x, Q):
    """All K roots, from the rational norm and both norm signs."""
    if x == 0:
        return [Q(0)]
    if x.b == 0:
        out = []
        for r, root in ((rational_sqrt(x.a), lambda t: Q(t)),
                        (rational_sqrt(x.a/5), lambda t: Q(0, t))):
            if r is not None:
                out.extend([root(r), -root(r)])
        return out
    norm = rational_sqrt(x.a*x.a-5*x.b*x.b)
    if norm is None:
        return []
    for n in (norm, -norm):
        u = rational_sqrt((x.a+n)/2)
        if u:
            root = Q(u, x.b/(2*u))
            if root*root == x:
                return [root, -root]
    return []


def check():
    sources = {p: digest(ROOT/p) for p in PINS}
    assert sources == PINS, 'pinned input changed'
    for p in (Path(__file__).relative_to(ROOT).as_posix(),
              'scripts/run_bounded_research_check_20261018.py'):
        sources[p] = digest(ROOT/p)
    Q = field(5)
    phi = Q(F(1, 2), F(1, 2))
    powers = [Q(1)]
    for _ in range(8):
        powers.append(powers[-1]*phi)
    zero = Q(0)
    unit = [Q(1)]+[zero]*6
    def dot(x, y):
        return -phi*x[0]*y[0]+sum((x[i]*y[i] for i in range(1, 7)), zero)
    def add(x, y):
        return [a+b for a, b in zip(x, y)]
    def scale(t, x):
        return [t*a for a in x]
    def columns(vectors):
        return [list(row) for row in zip(*vectors)]
    def mm(a, b):
        return [[sum((x*y for x, y in zip(row, col)), zero)
                 for col in zip(*b)] for row in a]
    def mv(a, x):
        return [sum((v*w for v, w in zip(row, x)), zero) for row in a]
    def enc(x):
        return [v.encoded() for v in x]
    def determinant(matrix):
        rows, ans = [list(row) for row in matrix], Q(1)
        for i in range(len(rows)):
            pivot = rows[i][i]
            assert pivot != 0
            ans = ans*pivot
            for j in range(i+1, len(rows)):
                t = rows[j][i]/pivot
                rows[j] = [a-t*b for a, b in zip(rows[j], rows[i])]
        return ans

    beta = []
    for i in range(1, 7):
        v = [zero]*7
        v[i] = Q(-1)
        if i < 6:
            v[i+1] = Q(1)
        beta.append(v)
    u = [phi-1, phi]+[zero]*5
    w = [powers[2]]+[phi]*5+[zero]
    h = [2*powers[3], 2*powers[2]]+[powers[3]]*5
    assert dot(u, u) == dot(w, w) == 2 and dot(h, h) == 1
    for outside, contacts in ((u, {0: -phi}), (w, {4: -phi}),
                              (h, {0: Q(-1), 5: -powers[3]})):
        assert all(dot(outside, v) == contacts.get(i, zero)
                   for i, v in enumerate(beta))
    assert dot(u, w) == dot(u, h) == dot(w, h) == 0
    aH, zH = h[0], sum(h[1:], zero)
    assert zH*zH-6*phi*aH*aH == 1 and zH == 12*phi+7
    tH = phi*aH/zH
    assert (phi-6*tH*tH).sign() > 0
    # The coefficient-box enclosure is exact, including the delicate B=14 cut.
    assert (Q(F(59*59, 4))-12*aH*aH).sign() > 0
    assert 5*140*140-313*313 == 31
    assert (Q(F(81, 25))-2*phi).sign() > 0
    profiles = [(1, 0, 'A5+A1'), (2, 0, 'A5+A1'),
                (2, 1, 'A6'), (2, 5, 'A6'),
                (2, 2, 'D6'), (2, 4, 'D6'), (2, 3, 'E6'),
                (1, 1, 'B6'), (1, 5, 'B6')]
    strips, binary, accepted, integral = [], {}, [], []
    tests = 0
    for r in (1, 3, 4, 6, 7, 12):
        times = []
        for b in range(14):
            for a in range(-1, 10):
                tests += 1
                time = a+b*phi
                if (time.sign() <= 0 or (r*aH*aH-time*time).sign() < 0
                    or (phi*r/6-time.conjugate()*time.conjugate()).sign() < 0):
                    continue
                times.append(enc([time])[0])
                for z in square_roots(6*phi*time*time+r, Q):
                    if z.sign() > 0:
                        assert z*z-6*phi*time*time == r
                        binary.setdefault(r, []).append((time, z))
        strips.append({'R': r, 'times': times})
    for q, k, kind in profiles:
        r = 6*q-k*(6-k) if k else 6*q
        for time, z in binary.get(r, []):
            y = (z-k)/6
            spatial = [y+int(i < k) for i in range(6)]
            v = [time]+spatial
            if not all(x.integral() for x in v):
                continue
            assert dot(v, v) == q
            assert (tH-phi*time/z).sign() >= 0
            assert all(dot(v, beta[i]) == (-1 if k == i+1 else 0)
                       for i in range(5))
            record = {'q': q, 'k': k, 'kind': kind, 'root': enc(v),
                      'beta6_product': dot(v, beta[5]).encoded(),
                      'u_product': dot(v, u).encoded()}
            integral.append(record)
            if all(dot(v, old).sign() <= 0 for old in beta+[u]):
                accepted.append(record)
    assert len(accepted) == 2, 'unknown first-wall candidates: '+json.dumps(
        {'accepted': accepted, 'integral_candidates': integral}, sort_keys=True)
    by_profile = {(row['q'], row['k']): row for row in accepted}
    assert set(by_profile) == {(1, 1), (2, 2)}
    g = add(scale(Q(2), h), beta[0])
    assert by_profile[(1, 1)]['root'] == enc(h), 'selected simple root differs from h'
    assert by_profile[(2, 2)]['root'] == enc(g), 'unknown nonsimple root'
    assert g == [4*powers[3]]+[4*phi+3]*2+[4*phi+2]*4
    assert dot(g, g) == 2 and g[0] == 2*aH and sum(g[1:], zero) == 2*zH
    assert phi*g[0]/sum(g[1:], zero) == tH
    assert h[2].integral() and h[2]*h[2].conjugate() == -1
    corner = beta[:5]+[h]
    gram = [[dot(v, w0) for w0 in corner] for v in corner]
    minors = [determinant([row[:k] for row in gram[:k]]) for k in range(1, 7)]
    assert all(x.sign() > 0 for x in minors) and minors[-1] == 1
    endpoint = [3*phi-1]+[Q(2)]*6
    assert dot(endpoint, endpoint) == -Q(1)/powers[7]
    assert all(dot(endpoint, v) == 0 for v in corner)
    # Exact full seven-dimensional reversal, not an assumed diagram symmetry.
    basis = columns([u]+beta)
    images = columns([w]+list(reversed(beta[:5]))+[h])
    amap = mm(images, inverse(basis, Q))
    assert all(x.integral() for row in amap for x in row)
    identity = [[Q(int(i == j)) for j in range(7)] for i in range(7)]
    assert mm(amap, amap) == identity
    assert all(dot(mv(amap, x), mv(amap, y)) == dot(x, y)
               for x in [unit]+beta for y in [unit]+beta)
    assert mv(amap, u) == w and mv(amap, w) == u and mv(amap, h) == beta[5]
    assert mv(amap, unit) == scale(powers[4], endpoint)
    assert mv(amap, unit)[0].sign() > 0
    assert sum((amap[i][i] for i in range(7)), zero) == 1
    def reflected(v):
        return add(v, scale(-2*dot(v, h), h))
    assert dot(beta[0], reflected(beta[0])) == 0
    assert dot(beta[5], reflected(beta[5])) == 1-2*powers[6]
    assert dot(beta[0], reflected(beta[5])) == -2*powers[3]
    assert 2*powers[6]-1 == 17+8*Q(0, 1)
    assert all(digest(ROOT/p) == pin for p, pin in sources.items()), 'source changed'
    return {'status': 'EXACT_ACTUAL_M6_EDGE_AND_REVERSAL_PASS', 'owns_lock': False,
            'sources': sources, 'field': 'Q(sqrt5); entries [rational,sqrt5 coefficient]',
            'coefficient_box': {'A': [-1, 9], 'B': [0, 13]}, 'tests': tests,
            'profiles': [{'q': q, 'k': k, 'type': kind} for q, k, kind in profiles],
            'strips': strips, 'integral_candidates': integral, 'accepted': accepted,
            'selected_simple': by_profile[(1, 1)],
            'nonsimple_identity': {'root': by_profile[(2, 2)],
                                   'positive_coefficients': {'h': 2, 'beta1': 1},
                                   'same_collision_parameter': True},
            'saturated_corner_principal_minors': enc(minors), 'endpoint': enc(endpoint),
            'endpoint_norm': dot(endpoint, endpoint).encoded(),
            'reversal_matrix': [enc(row) for row in amap], 'reversal_order': 2,
            'reversal_trace': '1', 'ambient_reversal_trace': '4',
            'written_inputs': ['complete norm1/2 pool and native A5 profiles',
                'actual old B6/u and sealed whole H3 interface',
                'actual first-wall existence before finite upper',
                'full-germ chamber identification', 'compact P6 double lineage',
                'no compact H9 conclusion or global infinite period']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    for flag in ('write', 'verify', 'stdout'):
        mode.add_argument('--'+flag, action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error('refusing to overwrite existing certificate')
    signal.alarm(45)
    data = check()
    payload = json.dumps(data, sort_keys=True, indent=2)+'\n'
    assert len(payload.encode()) <= 131072
    if args.write:
        with CERT.open('x') as stream:
            stream.write(payload)
    elif args.verify:
        assert CERT.stat().st_size <= 131072
        assert data == json.loads(CERT.read_text()), 'saved certificate differs'
    print(payload if args.stdout else json.dumps({'status': 'PASS', 'profiles': 9,
          'accepted': 2, 'simple': 1, 'mode': 'write' if args.write else 'verify'}),
          end='' if args.stdout else '\n')


if __name__ == '__main__':
    main()
