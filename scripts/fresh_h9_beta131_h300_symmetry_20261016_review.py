#!/usr/bin/env python3
"""Independent two-corner proof of nonreflectivity of the norm131 host.

No producer/geometry/arithmetic-helper imports. Project Golden arithmetic,
literal form, direct orthogonal projections and reconstructed frame map.
"""
import argparse
import json
import sys
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from coxeter_polytopes.exact.golden import Golden as K

TAG = 'fresh_h9_beta131_h300_symmetry_20261016_review'
PREFIX = 'results/fresh_h9_beta131_h300_20261016.json'
CENSUS_REVIEW = 'results/fresh_h9_beta131_h300_20261016_review.json'
PRIMARY = 'results/fresh_h9_beta131_h300_20261016_symmetry.json'
LEFT = [0, 1, 2, 3, 4, 5, 6, 7, 10]
RIGHT = [5, 2, 4, 3, 1, 244, 42, 11, 7]


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def read(path):
    obj = json.loads((ROOT / path).read_text())
    for source, digest in obj.get('sources', {}).items():
        assert sha(source) == digest, (path, source)
    return obj


def enc(x):
    if isinstance(x, K):
        return x.to_json()
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [enc(v) for v in x]
    return x


def dec(v):
    return [K.from_json(x) for x in v]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), K())


def mv(a, x):
    return [dot(row, x) for row in a]


def trans(a):
    return list(map(list, zip(*a)))


def mm(a, b):
    return [[dot(row, col) for col in zip(*b)] for row in a]


def identity(n):
    return [[K(int(i == j)) for j in range(n)] for i in range(n)]


def inverse(a):
    n = len(a)
    rows = [list(row)+identity(n)[i] for i, row in enumerate(a)]
    determinant = K(1)
    for j in range(n):
        i = next(i for i in range(j, n) if rows[i][j])
        if i != j:
            rows[i], rows[j] = rows[j], rows[i]
            determinant = -determinant
        d = rows[j][j]
        determinant *= d
        rows[j] = [x/d for x in rows[j]]
        for i in range(n):
            if i != j:
                d = rows[i][j]
                rows[i] = [x-d*y for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows], determinant


def pivots(a):
    a = [row[:] for row in a]
    out = []
    while a:
        pivot = a[0][0]
        assert pivot > 0
        out.append(pivot)
        a = [[a[i][j]-a[i][0]*a[0][j]/pivot
              for j in range(1, len(a))] for i in range(1, len(a))]
    return out


def run(require_prefix_review=True):
    prefix, saved = read(PREFIX), read(PRIMARY)
    assert prefix['height_bound'] == '300' and len(prefix['roots']) == 282
    phi = K(0, 1)
    g = [[K() for j in range(10)] for i in range(10)]
    g[0][0], g[0][9], g[9][0], g[9][9] = -2*phi, K(-1), K(-1), K(6)
    edges = {(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)}
    for i in range(8):
        for j in range(8):
            g[i+1][j+1] = K(2 if i == j else -1 if (min(i, j), max(i, j)) in edges else 0)
    assert g == [dec(row) for row in prefix['form']]
    _, determinant = inverse(g)
    assert determinant == -(1+12*phi)
    assert determinant.norm() == -131
    conjugate_positive = pivots([[x.conjugate() for x in row] for row in g])
    assert len(conjugate_positive) == 10
    roots = [dec(row) for row in prefix['roots']]
    base = [K(6)]+[K()]*8+[K(1)]
    assert base == dec(prefix['basepoint'])

    def pair(x, y):
        return dot(x, mv(g, y))

    assert pair(base, base) == -6*(1+12*phi)
    base_pairs = [pair(base, n) for n in roots]
    assert all(x <= 0 for x in base_pairs)
    assert [i for i, x in enumerate(base_pairs) if not x] == list(range(8))
    center = dec(prefix['strict_center'])
    assert pair(center, center) < 0 and pair(center, base) < 0
    assert all(pair(center, n) < 0 for n in roots)
    corners = []
    common = None
    for ids, key in ((LEFT, 'source_corner'), (RIGHT, 'target_corner')):
        ns = [roots[i] for i in ids]
        s = [[pair(a, b) for b in ns] for a in ns]
        spherical_pivots = pivots(s)
        assert all(s[i][j] <= 0 for i in range(9) for j in range(i))
        j, dets = inverse(s)
        assert all(x >= 0 for row in j for x in row)
        coefficients = mv(j, [base_pairs[i] for i in ids])
        assert all(x <= 0 for x in coefficients)
        # Compute each apex from the basepoint, independently of saved points.
        v = [base[k]-sum((coefficients[i]*ns[i][k] for i in range(9)), K()) for k in range(10)]
        assert all(pair(v, n) == 0 for n in ns)
        norm, future = pair(v, v), pair(v, base)
        assert norm < 0 and future == norm
        pairings = [pair(v, n) for n in roots]
        assert all(x <= 0 for x in pairings)
        assert [i for i, x in enumerate(pairings) if not x] == sorted(ids)
        raw = dec(saved[key]['point'])
        scale = pair(base, raw)/pair(raw, raw)
        assert scale > 0 and v == [scale*x for x in raw]
        if common is None:
            common = s
        else:
            assert s == common
        corners.append(dict(ids=ids, apex=v, norm=norm, coefficients=coefficients,
                            pairings=pairings, spherical_pivots=spherical_pivots,
                            comparison_to_saved_apex=scale))
    v, w = [c['apex'] for c in corners]
    frame = trans([roots[i] for i in LEFT]+[v])
    other = trans([roots[i] for i in RIGHT]+[w])
    fi, ds = inverse(frame)
    _, dt = inverse(other)
    scale = ds/dt
    if scale < 0:
        scale = -scale
    assert scale > 0 and scale*scale*pair(w, w) == pair(v, v)
    for row in other:
        row[-1] *= scale
    matrix = mm(other, fi)
    inv, det = inverse(matrix)
    assert det == -1
    assert all(x.is_integral for row in matrix+inv for x in row)
    assert mm(trans(matrix), mm(g, matrix)) == g
    assert mm(matrix, inv) == mm(inv, matrix) == identity(10)
    assert matrix == [dec(row) for row in saved['matrix']]
    assert inv == [dec(row) for row in saved['inverse_matrix']]
    for i, j in zip(LEFT, RIGHT):
        assert mv(matrix, roots[i]) == roots[j]
    assert mv(matrix, v) == [scale*x for x in w]
    assert pair(v, w) < 0  # positive apex scaling preserves the future sheet.
    trace = sum((matrix[i][i] for i in range(10)), K())
    assert trace == 32+53*phi and trace > 10
    paths = [f'scripts/{TAG}.py', 'src/coxeter_polytopes/exact/golden.py', PREFIX, PRIMARY,
             'results/fresh_h9_beta131_h300_20261016.md']
    if require_prefix_review:
        review = read(CENSUS_REVIEW)
        assert review['controls'] == 200 and review['walls'] == 282
        assert review['height_bound'] == '300'
        assert review['status'] == 'PASS_INDEPENDENT_COMPLETE_BETA131_PREFIX'
        paths.append(CENSUS_REVIEW)
    return enc(dict(status='PASS_INDEPENDENT_BETA131_FULL_CHAMBER_NONREFLECTIVITY' if require_prefix_review
                    else 'PASS_LOCAL_SYMMETRY_PENDING_INDEPENDENT_PREFIX_REVIEW',
                    sources={p: sha(p) for p in paths}, corners=corners,
                    form=g, determinant=determinant, conjugate_positive_pivots=conjugate_positive,
                    basepoint=base, common_spherical_gram=common, positive_apex_scale=scale,
                    matrix=matrix, inverse_matrix=inv, matrix_determinant=det, trace=trace,
                    conclusion='The norm131 lattice has an infinite-order integral future symmetry of its full reflection chamber and is nonreflective.',
                    full_conclusion_verified=require_prefix_review,
                    scope='Only two literal finite corner frames are needed. No reliance on completeness of the 70-corner census or the corner-map search. Completeness of the independent bounded Vinberg prefix and all primitive norm classes is required for the full-chamber conclusion.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--witness-only', action='store_true')
    args = ap.parse_args()
    assert not (args.verify and args.witness_only)
    result = run(not args.witness_only)
    if not args.witness_only:
        path = ROOT / 'results' / f'{TAG}.json'
        if args.verify:
            assert result == json.loads(path.read_text())
        else:
            path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['trace'])
