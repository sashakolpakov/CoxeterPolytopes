#!/usr/bin/env python3
"""Independent exact geometry check; imports no producer arithmetic.

Each normal is A+B*sqrt(D), over rational pairs Q(sqrt(5)). A mutual
pairing has two extra positive radicals, handled by nested exact signs.
"""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_makarov_joint_20261013_exact_independent'
SOURCE = ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
SEED = ROOT/'results/fresh_h9_makarov_face_extension_20261013.json'


def k(a=0, b=0):
    return F(a), F(b)


Z = k()


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    a, b = x
    c, d = y
    return a*c+5*b*d, a*d+b*c


def sign(x):
    a, b = x
    if not b:
        return (a > 0)-(a < 0)
    if not a:
        return (b > 0)-(b < 0)
    if (a > 0) == (b > 0):
        return (a > 0)-(a < 0)
    d = a*a-5*b*b
    return ((a > 0)-(a < 0))*((d > 0)-(d < 0))


def reciprocal(x):
    a, b = x
    d = a*a-5*b*b
    assert d
    return a/d, -b/d


def inertia(Q):
    M = [row[:] for row in Q]
    positive = negative = 0
    while M:
        pivot = M[0][0]
        sg = sign(pivot)
        assert sg
        positive += sg > 0
        negative += sg < 0
        M = [[sub(M[i][j], mul(mul(M[i][0], M[0][j]), reciprocal(pivot)))
              for j in range(1, len(M))] for i in range(1, len(M))]
    return [positive, negative]


def lmul(x, y, D):
    a, b = x
    c, d = y
    return add(mul(a, c), mul(D, mul(b, d))), add(mul(a, d), mul(b, c))


def lsign(x, D):
    a, b = x
    sa, sb = sign(a), sign(b)
    if not sign(D) or not sb:
        return sa
    assert sign(D) > 0
    if not sa:
        return sb
    if sa == sb:
        return sa
    return sa*sign(sub(mul(a, a), mul(D, mul(b, b))))


def psign(x, D, E):
    a, b = x
    sa, sb = lsign(a, D), lsign(b, D)
    if not sign(E) or not sb:
        return sa
    assert sign(E) > 0
    if not sa:
        return sb
    if sa == sb:
        return sa
    aa, bb = lmul(a, a, D), lmul(b, b, D)
    return sa*lsign(tuple(sub(aa[i], mul(E, bb[i])) for i in range(2)), D)


def pmul(x, y, D, E):
    ac, bd = lmul(x[0], y[0], D), lmul(x[1], y[1], D)
    ad, bc = lmul(x[0], y[1], D), lmul(x[1], y[0], D)
    return (tuple(add(ac[i], mul(E, bd[i])) for i in range(2)),
            tuple(add(ad[i], bc[i]) for i in range(2)))


def scalar(x):
    return ((x, Z), (Z, Z))


def psub(x, y):
    return tuple(tuple(sub(x[i][j], y[i][j]) for j in range(2)) for i in range(2))


def pair(a, b, Q):
    # Normal a uses the first radical, b the second radical.
    out = [[Z, Z], [Z, Z]]
    for i in range(len(a)):
        for j in range(len(b)):
            for s in range(2):
                for t in range(2):
                    out[t][s] = add(out[t][s], mul(Q[i][j], mul(a[i][s], b[j][t])))
    return tuple(tuple(row) for row in out)


def encode(x):
    return [[list(map(str, a)) for a in row] for row in x]


def run(verify=False):
    source = json.loads(SOURCE.read_text())
    seed = json.loads(SEED.read_text())
    for key in ['form', 'seed_normals', 'strict_center']:
        assert source[key] == seed[key]
    Q = [[k(*a) for a in row] for row in seed['form']]
    assert inertia(Q) == [9, 1]
    normals = [[(k(*a), Z) for a in row] for row in seed['seed_normals']]
    D = [Z]*len(normals)
    caps = source['caps']
    for row in caps:
        d = k(*row['radicand'])
        assert sign(d) > 0
        D.append(d)
        normals.append([(k(*a), k(*b)) for a, b in row['normal']])
    n = len(normals)
    expected_norm = [k(2)]*n
    expected_norm[6] = k(F(7, 11), F(-1, 11))
    for i in range(n):
        p = pair(normals[i], normals[i], Q)
        assert psign(psub(p, scalar(expected_norm[i])), D[i], D[i]) == 0
    center = [(k(*a), Z) for a in seed['strict_center']]
    assert psign(pair(center, center, Q), Z, Z) < 0
    center_pairings = []
    for i, r in enumerate(normals):
        p = pair(r, center, Q)
        assert psign(p, D[i], Z) < 0
        center_pairings.append(encode(p))

    alphabet = {3: k(F(1, 4)), 4: k(F(1, 2)),
                5: k(F(3, 8), F(1, 8)), 6: k(F(3, 4)),
                10: k(F(5, 8), F(1, 8))}
    pairs = []
    counts = Counter()
    for i in range(n):
        for j in range(i):
            p = pair(normals[i], normals[j], Q)
            sg = psign(p, D[i], D[j])
            assert sg <= 0, (i, j, 'obtuse')
            if not sg:
                tag = 'm2'
            else:
                square = pmul(p, p, D[i], D[j])
                product = mul(expected_norm[i], expected_norm[j])
                gap = psign(psub(square, scalar(product)), D[i], D[j])
                if gap >= 0:
                    tag = 'parallel' if not gap else 'dotted'
                else:
                    labels = [m for m, c in alphabet.items()
                              if psign(psub(square, scalar(mul(product, c))), D[i], D[j]) == 0]
                    assert len(labels) == 1, (i, j, 'unrecognized exact angle')
                    tag = 'm'+str(labels[0])
            counts[tag] += 1
            pairs.append({'walls': [j, i], 'type': tag, 'pairing': encode(p)})

    escape_int = [140532998, 178861908, 130353897, 81870848, 33432956,
                  12307295, -284558239, -345319008, -230212893, -115106557]
    escape = [(k(a), Z) for a in escape_int]
    escape_norm = pair(escape, escape, Q)
    assert psign(escape_norm, Z, Z) > 0
    assert psign(pair(escape, center, Q), Z, Z) < 0
    escape_pairings = []
    for i, r in enumerate(normals):
        p = pair(r, escape, Q)
        assert psign(p, D[i], Z) < 0
        escape_pairings.append(encode(p))

    # Independently verify the retained old corners and all forty first
    # finite endpoints, including their full wall signs and actual ranks.
    old_vertices = [[(k(*a), Z) for a in row] for row in seed['compact_H5_vertices']]
    old_vertices = [row+[(Z, Z)]*4 for row in old_vertices]
    for vertex, active in zip(old_vertices, seed['vertex_active_sets'], strict=True):
        assert psign(pair(vertex, vertex, Q), Z, Z) < 0
        actual = []
        for i, r in enumerate(normals):
            sg = psign(pair(vertex, r, Q), Z, D[i])
            assert sg <= 0
            if not sg:
                actual.append(i)
        assert actual == active+list(range(7, 11))
        gram = [[pair(normals[i], normals[j], Q)[0][0] for j in actual] for i in actual]
        assert inertia(gram) == [9, 0]
    B = [[k(*a) for a in row] for row in seed['H4_inverse']]
    seen = set()
    endpoint_records = []
    for row in source['new_finite_endpoints']:
        d = k(*row['radicand'])
        assert sign(d) > 0
        point = [(k(*a), k(*b)) for a, b in row['point']]
        parameter = tuple(k(*a) for a in row['edge_parameter'])
        assert lsign(parameter, d) > 0
        vi, released, owner = row['old_vertex'], row['released_H4_wall'], row['new_owner']
        assert 0 <= vi < 10 and 7 <= released < 11 and 11 <= owner < n
        vertex = old_vertices[vi]
        for i in range(6):
            assert lsign(tuple(sub(point[i][a], vertex[i][a]) for a in range(2)), d) == 0
        for i in range(4):
            predicted = tuple(neg(mul(a, B[i][released-7])) for a in parameter)
            assert lsign(tuple(sub(point[i+6][a], predicted[a]) for a in range(2)), d) == 0
        assert psign(pair(point, point, Q), d, d) < 0
        assert psign(pair(point, center, Q), d, Z) < 0
        incident = []
        for i, r in enumerate(normals):
            sg = psign(pair(point, r, Q), d, D[i])
            assert sg <= 0
            if not sg:
                incident.append(i)
        support = seed['vertex_active_sets'][vi]+[i for i in range(7, 11) if i != released]
        assert support == row['support']
        assert incident == sorted(support+[owner]) == row['incident']
        gram = [[pair(normals[i], normals[j], Q)[0][0] for j in support] for i in support]
        assert inertia(gram) == [8, 0]
        # Strictness at the old vertex proves the owner is independent of
        # its eight old incident walls. Both endpoint signs then prove
        # firstness along the segment by linearity of every inequality.
        assert psign(pair(vertex, normals[owner], Q), Z, D[owner]) < 0
        seen.add(tuple(incident))
        endpoint_records.append({'old_vertex': vi, 'released_H4_wall': released,
                                 'incident': incident, 'actual_rank': 9})
    assert len(endpoint_records) == len(seen) == 40
    assert {(r['old_vertex'], r['released_H4_wall']) for r in endpoint_records} == {
        (vi, wall) for vi in range(10) for wall in range(7, 11)}
    out = {'status': 'PASS_INDEPENDENT_EXACT_JOINT_CHAMBER_PAIRS_AND_ESCAPE',
           'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
           'seed_sha256': hashlib.sha256(SEED.read_bytes()).hexdigest(),
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'arithmetic': 'Standalone Fraction rational pairs and two nested quadratic sign comparisons; no producer imports.',
           'facet_count': n, 'pair_count': len(pairs), 'pair_types': dict(sorted(counts.items())),
           'independent_form_inertia': [9, 1],
           'all_norms_verified': True, 'strict_center_pairings': center_pairings,
           'all_pairs': pairs, 'strict_integer_escape': escape_int,
           'escape_norm': encode(escape_norm), 'strict_escape_pairings': escape_pairings,
           'retained_old_finite_vertices': 10,
           'distinct_first_finite_endpoints': endpoint_records,
           'actual_facets': 'For each r, c-(c.r)/(r.r)*r is timelike, lies on r, and pairs strictly negatively with every other wall by mutual acuteness.',
           'scope': 'Verifies form signature, all exact wall norms and Coxeter pairings, a common strict timelike center, actual facets, a strict spacelike escape, ten retained old finite vertices, and forty distinct actual first finite endpoints. Does not independently enumerate the search recipes or all new vertices. No finite-volume example.'}
    path = ROOT/'results'/f'{TAG}.json'
    if verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], n, 'facets;', len(pairs), 'pairs;', dict(counts))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    run(ap.parse_args().verify)
