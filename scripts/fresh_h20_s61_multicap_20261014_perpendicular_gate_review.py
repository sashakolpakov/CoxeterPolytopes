#!/usr/bin/env python3
"""Independent small arithmetic/proof audit of the perpendicular S61 gate.

Uses the earlier independent sparse field and dyadic interval implementation;
imports no primary construction or gate arithmetic.  Does not repeat the
already independently certified 133-wall pair and polygon census.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field, parse
import fresh_h20_bridge_faces_20261014_review as B
import fresh_h20_s61_multicap_20261014_exact_review as I

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h20_s61_multicap_20261014_perpendicular_gate'
TAG = BASE + '_review'
FILES = [
    'results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json',
    'results/fresh_h20_e7_bridge_face_20261014_exact.json',
    'results/fresh_h20_bridge_edge_closure_20261014.json',
    'results/fresh_h20_s61_multicap_20261014_exact.json',
    'results/fresh_h20_s61_multicap_20261014_exact_review.json',
    f'results/{BASE}.json',
]


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def strings(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, (tuple, list)):
        return [strings(y) for y in x]
    if isinstance(x, dict):
        return {k: strings(v) for k, v in x.items()}
    return x


def run():
    old, bridge, edge, full, previous, gate = [
        json.loads((ROOT / p).read_text()) for p in FILES]
    for doc in [old, bridge, edge, full, previous, gate]:
        for p, h in doc['sources'].items():
            assert sha(p) == h, p
    assert previous['status'] == 'PASS_INDEPENDENT_133_WALL_RIGHT_CHAIN_AND_COMPACT_S61_FACE'
    assert previous['polygon_vertices'] == 29
    assert previous['nonconsecutive_pair_exclusions'] == 377
    ff = Field([parse(old['field']['discriminant'], 0),
                parse(bridge['field']['new_positive_sqrt_radical'], 1)])
    add, mul, sc, total, c = ff.add, ff.mul, ff.scale, ff.total, ff.const
    zero = ff.zeros[2]
    Q = [list(map(F, row)) for row in old['form']]
    assert B.inertia(Q) == [20, 1, 0]
    rows = [[(j, q) for j, q in enumerate(row) if q] for row in Q]
    cache = {}

    def vp(raw, level):
        return tuple(ff.lift(parse(x, level), level) for x in raw)

    def dot(v, w):
        if id(w) not in cache:
            cache[id(w)] = (w, tuple(total(sc(w[j], q) for j, q in row)
                                      for row in rows))
        return total(mul(x, y) for x, y in zip(v, cache[id(w)][1]))

    def rational(x):
        for level in [2, 1]:
            assert x[1] == ff.zeros[level - 1]
            x = x[0]
        assert not any(x[1:])
        return x[0]

    @lru_cache(None)
    def bounds(x):
        return I.out(ff.interval(x, 2, I.BITS))

    N = [vp(v, 1) for v in old['normals']] + [vp(bridge['new_normal'], 2)]
    S = full['support']
    assert S == edge['recipe']['support'] == gate['support']
    assert len(S) == 18 and edge['recipe']['contact_wall'] == 61
    roots = [N[i] for i in S]
    G = [[rational(dot(v, w)) for w in roots] for v in roots]
    assert G == [list(map(F, row)) for row in full['support_gram']]
    B.pd(G)
    assert all(x >= 0 for row in B.inv(G) for x in row)
    X, r, W = [vp(full['frame'][k], 2) for k in ['X', 'r95', 'W']]
    assert X == vp(edge['recipe']['X'], 2)
    assert r == N[95] and W == vp(edge['recipe']['W'], 2)
    A, w = [parse(full['frame'][k], 2) for k in ['A', 'W_norm']]
    assert [[dot(v, u) for u in [X, r, W]] for v in [X, r, W]] == [
        [sc(A, -1), zero, zero], [zero, c(2), zero], [zero, zero, w]]
    assert ff.sign(A) > 0 and ff.sign(add(w, c(-2))) > 0
    assert all(dot(v, n) == zero for v in [X, r, W] for n in roots)
    E = tuple(c(x) for x in full['strict_escape']['rational_point'])
    qE = rational(dot(E, E))
    assert qE == F(38017800175279, 500000000000)
    assert full['strict_escape']['rational_point'] == gate['original_escape']
    den = [I.root(bounds(x)) for x in [A, c(2), w]]
    covector = [I.div(bounds(dot(v, E)), d) for v, d in zip([X, r, W], den)]
    projection = [I.neg(covector[0]), covector[1], covector[2]]
    qproj = I.lor(projection, projection)
    assert projection[0][0] > 0 and qproj[1] < -76
    a = tuple(map(F, full['parameter']['isolating_interval']))
    coeff, _ = I.chain(a)
    rw, r2 = den[2], den[1]
    k = I.root(I.minus(bounds(w), I.pt(F(5, 6))))
    sides = [[I.pt(0), r2, I.pt(0)],
             [k, I.neg(I.div(I.pt(1), r2)), rw]] + coeff + [
             [k, I.neg(I.div(I.pt(1), r2)), I.neg(rw)]]
    ids = [95, 96] + list(range(107, 133)) + [97]
    products = [I.lor(projection, v) for v in sides]
    assert len(products) == 29 and all(x[1] < -15 for x in products)
    for actual, saved in zip(projection, gate['projected_escape_frame_bounds']):
        assert I.overlap(actual, tuple(map(F, saved)))
    assert I.overlap(qproj, tuple(map(F, gate['projected_escape_norm_bound'])))
    for actual, saved in zip(products, gate['projected_escape_polygon_pair_bounds']):
        assert I.overlap(actual, tuple(map(F, saved)))

    # Rebuild all actual PD20 frames, with each appended root normalized to q=2.
    # The full polygon census and the physical normalization are already bound
    # by the separate preceding review; these are fresh rational matrix checks.
    expected_pairs = {tuple(sorted([ids[j], ids[(j + 1) % 29]]))
                      for j in range(29)}
    assert expected_pairs == {tuple(sorted(row['walls'])) for row in full['face']['vertices']}
    classes = Counter()
    checks = []
    for first, second in sorted(expected_pairs):
        H = [row + [F(0), F(0)] for row in G] + [[F(0)] * 20 for _ in range(2)]
        for j, wall in enumerate([first, second], 18):
            H[j][j] = 2
            if wall in [96, 97]:
                H[j][S.index(61)] = H[S.index(61)][j] = -1
        H[18][19] = H[19][18] = -int(95 in [first, second])
        pivots = B.pd(H)
        J = B.inv(H)
        assert all(x >= 0 for row in J for x in row)
        kind = ('old95_corner' if 95 in [first, second] else
                'chain_endpoint' if {first, second} & {96, 97} else 'chain_interior')
        classes[kind] += 1
        checks.append(dict(walls=[first, second], minimum_inverse=min(x for row in J for x in row),
                           minimum_positive_pivot=min(pivots)))
    assert dict(classes) == {'old95_corner': 2, 'chain_endpoint': 2, 'chain_interior': 25}
    paths = FILES + [f'scripts/{BASE}.py', f'scripts/{TAG}.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact_review.py',
        'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
        'scripts/fresh_h20_bridge_faces_20261014_review.py']
    return strings(dict(
        status='PASS_INDEPENDENT_ARITHMETIC_AND_PROOF_PERPENDICULAR_S61_GATE',
        sources={p: sha(p) for p in paths}, interval_bits=I.BITS, support=S,
        original_escape_norm=qE, projected_escape_frame_bounds=projection,
        projected_escape_norm_bound=qproj, all_29_polygon_pair_bounds=products,
        active_frame_counts=dict(classes), actual_vertex_Gram_checks=checks,
        proof_checks=dict(future_projection_sign=True, common_future_strict_center_required=True,
                          inverse_positive_vertex_preservation=True,
                          compact_Klein_convex_hull=True, finite_extension_infinite_volume=True),
        scope='All real spacelike added normals perpendicular to S61 and acute with the retained 133 walls, with a common future timelike strict interior. No label/field/norm cutoff. A finite family leaves the strict spacelike escape; nonright S61 contacts and changed cores are outside this theorem.',
        independence='No producer arithmetic imports. Independently reconstructs the physical frame, escape projection, 29 strict side signs and 29 rational PD20 inverse checks. Inherits the separately independently replayed complete compact polygon and old 133-wall geometry; no repeated full pair/corner census.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], flush=True)
