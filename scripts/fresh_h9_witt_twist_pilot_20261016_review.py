#!/usr/bin/env python3
"""Independent meet-in-the-middle census for the locally distinct H9 pilot.

Enumerates complete ordered spatial norm shells before chronological wall
filtering, unlike the producer's pruned coordinate recursion.
"""
import argparse
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
from math import ceil, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from coxeter_polytopes.exact.golden import Golden as K

TAG = 'fresh_h9_witt_twist_pilot_20261016_review'
PRIMARY = 'results/fresh_h9_witt_twist_pilot_20261016.json'
SETUP_REVIEW = 'results/fresh_h9_witt_twist_20261016_review.json'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def enc(x):
    if isinstance(x, K):
        return x.to_json()
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    return x


def dec(v):
    return [K.from_json(x) for x in v]


def key(v):
    return tuple((x.a, x.b) for x in v)


def sqrtceil(x):
    a = isqrt(x.numerator // x.denominator)
    return a+int(a*a < x)


def coefficient_pool(physical, conjugate):
    # Deliberately larger boxes than the producer; exact filters follow.
    aa = ceil(Q(17*conjugate+7*physical, 22))+2
    bb = ceil(Q(5*(physical+conjugate), 11))+2
    return [K(a, b) for a in range(-aa, aa+1) for b in range(-bb, bb+1)], [aa, bb]


def run():
    data, setup = [json.loads((ROOT / path).read_text()) for path in (PRIMARY, SETUP_REVIEW)]
    for document in (data, setup):
        for path, digest in document['sources'].items():
            assert sha(path) == digest, path
    assert setup['status'] == 'PASS_INDEPENDENT_DISTINCT_WITT_HOST_NORMS_AND_ZERO_SECTOR'
    cap = Q(data['height_bound'])
    phi, beta = K(0, 1), K(4, -1)
    diagonal = [-phi]+[K(1)]*7+[beta]*2
    assert [dec(row) for row in data['form']] == [[diagonal[i] if i == j else K() for j in range(10)] for i in range(10)]
    norms = [K(1), K(2), beta, 2*beta]
    assert norms == dec(data['norm_classes']) == dec(setup['primitive_norm_classes'])
    assert phi > K(Q(8, 5)) and -phi.conjugate() > K(Q(3, 5))
    assert phi < K(Q(17, 10)) and phi-1 < K(Q(7, 10)) and 2*phi-1 > K(Q(11, 5))
    assert beta > 2 and beta.conjugate() > 4
    assert all(q < 5 and q.conjugate() < 10 for q in norms)
    xp, xb = coefficient_pool(sqrtceil(Q(25, 8)*cap), 5)
    yp, yb = coefficient_pool(sqrtceil(5*(1+cap)), 4)
    zp, zb = coefficient_pool(sqrtceil(Q(5, 2)*(1+cap)), 2)
    controls = []
    for q in norms:
        for x in xp:
            if x <= 0 or (q in norms[2:] and not (x/beta).is_integral):
                continue
            h, target = phi*x*x/q, q+phi*x*x
            if h <= cap and target.conjugate() >= 0:
                controls.append((h, q, x))
    controls.sort(key=lambda row: row[0])
    assert len({h for h, q, x in controls}) == len(controls)
    assert controls == [(K.from_json(r['height']), K.from_json(r['norm']),
                         K.from_json(r['control'])) for r in data['levels']]

    def pair(a, b):
        return sum((d*x*y for d, x, y in zip(diagonal, a, b)), K())

    normals = [dec(row) for row in setup['zero_simple_walls']]
    assert normals == [dec(row) for row in data['zero_simple_walls']]
    records, totals = [], Counter()
    for index, (height, q, x) in enumerate(controls):
        target = q+phi*x*x
        def admissible(value, weight, is_y):
            if value < 0 or (is_y and q in norms[2:] and not (value/beta).is_integral):
                return False
            remainder = target-weight*value*value
            return remainder >= 0 and remainder.conjugate() >= 0
        ys = sorted([v for v in yp if admissible(v, K(1), True)], reverse=True)
        zs = sorted([v for v in zp if admissible(v, beta, False)], reverse=True)
        saved = data['levels'][index]
        assert set(ys) == set(dec(saved['y_alphabet'])) and set(zs) == set(dec(saved['z_alphabet']))
        counts = Counter()
        def block(values, length, weight):
            squares = {z: weight*z*z for z in values}
            out = []
            for v in combinations_with_replacement(values, length):
                counts[f'block{length}_candidates'] += 1
                energy = sum((squares[z] for z in v), K())
                if energy <= target and energy.conjugate() <= target.conjugate():
                    out.append((v, energy))
            return out
        first, second, zblock = block(ys, 3, K(1)), block(ys, 4, K(1)), block(zs, 2, beta)
        by_energy = defaultdict(list)
        for v, energy in second:
            by_energy[energy].append(v)
        found = []
        for z, ez in zblock:
            for left, el in first:
                remainder = target-ez-el
                for right in by_energy.get(remainder, ()):
                    if left[-1] < right[0]:
                        continue
                    counts['ordered_norm_shell'] += 1
                    r = [x, *left, *right, *z]
                    assert all(v.is_integral for v in r) and pair(r, r) == q
                    assert all((2*d*v/q).is_integral for d, v in zip(diagonal, r))
                    if all(pair(r, n) <= 0 for n in normals):
                        found.append(r)
        assert len({key(r) for r in found}) == len(found)
        assert {key(r) for r in found} == {key(dec(r)) for r in saved['roots']}
        assert all(pair(r, s) <= 0 for i, r in enumerate(found) for s in found[:i])
        normals.extend(found)
        totals.update(counts)
        records.append(dict(height=height, norm=q, control=x, counts=dict(counts),
                            roots=sorted(found, key=key)))
        print('independent Witt pilot', index+1, '/', len(controls), 'walls', len(normals), flush=True)
    assert {key(r) for r in normals} == {key(dec(r)) for r in data['roots']}
    normals = [dec(r) for r in data['roots']]
    gram = [[pair(r, s) for s in normals] for r in normals]
    assert gram == [dec(row) for row in data['gram']]
    assert all(gram[i][j] <= 0 for i in range(len(normals)) for j in range(i))
    h = [K(1)]+[K()]*9
    center = dec(data['strict_center'])
    assert pair(center, center) < 0 and pair(center, h) < 0
    assert all(pair(center, r) < 0 for r in normals)
    assert len(data['actual_facet_points']) == len(normals)
    for i, raw in enumerate(data['actual_facet_points']):
        p = dec(raw)
        assert pair(p, p) < 0 and pair(p, h) < 0
        assert all(pair(p, r).sign() == (0 if j == i else -1) for j, r in enumerate(normals))
    pairs = Counter()
    finite = {K(): 'right', K(1): 'pi/3', K(2): 'pi/4', 1+phi: 'pi/5',
              K(3): 'pi/6', 2+phi: 'pi/10'}
    for i in range(len(normals)):
        for j in range(i):
            cycle = 4*gram[i][j]*gram[i][j]/(gram[i][i]*gram[j][j])
            assert cycle != 4
            pairs['dotted' if cycle > 4 else finite[cycle]] += 1
    assert dict(pairs) == data['pair_counts']
    paths = [f'scripts/{TAG}.py', 'src/coxeter_polytopes/exact/golden.py', PRIMARY, SETUP_REVIEW]
    return enc(dict(status='PASS_INDEPENDENT_COMPLETE_WITT_TWIST_PREFIX',
                    sources={p: sha(p) for p in paths}, height_bound=cap,
                    independent_coefficient_boxes=dict(control=xb, y=yb, z=zb),
                    controls=len(controls), walls=len(normals), counts=dict(totals), levels=records,
                    pair_counts=dict(pairs),
                    scope='Complete bounded root census by independent spatial-shell joins, plus every facet and Gram witness. No complete cone census or finite-volume/reflectivity assertion.'))


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
    print(result['status'], result['controls'], result['walls'], result['counts'])
