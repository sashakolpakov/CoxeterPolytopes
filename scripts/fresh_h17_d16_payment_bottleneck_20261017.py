#!/usr/bin/env python3
"""Exact local-payment demands in the frozen finite-volume H17 pyramid.

No new geometry, dimension bound, or optimization.  Reconstruct all 136
two-faces at its simple D16 cusp from the saved complete chamber cone.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h17_d16_payment_bottleneck_20261017'
OLD = 'fresh_nonincident_carrier_incidence_20261017'
MIXED = 'fresh_ideal_distance7_mixed_replacement_20261017'
NODES = set(range(19))
LEFT, RIGHT = set(range(9)), set(range(10, 19))
EDGES = [(0, 1), (0, 2), (0, 4), (2, 3), (4, 5),
         (5, 6), (6, 7), (7, 8), (8, 9)]
EDGES += [(a + 10 if a < 9 else 9, b + 10 if b < 9 else 9)
          for a, b in EDGES.copy()]
GRAPH = {i: set() for i in NODES}
for a, b in EDGES:
    GRAPH[a].add(b)
    GRAPH[b].add(a)
CUSP = NODES - {3, 13}


def sha(p):
    h = sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def component_and_distances(nodes, start):
    nodes = set(nodes)
    d, q = {start: 0}, deque([start])
    while q:
        a = q.popleft()
        for b in sorted((GRAPH[a] & nodes) - d.keys()):
            d[b] = d[a] + 1
            q.append(b)
    return set(d), d


def histogram(rows, key):
    return {str(k): v for k, v in sorted(Counter(key(r) for r in rows).items())}


def run():
    old = json.loads((ROOT / f'results/{OLD}.json').read_text())
    vertices = old['vertices']
    assert len(vertices) == 82
    assert sum(v['ideal'] for v in vertices) == 2
    assert [v['support'] for v in vertices if v['ideal']] == [[9], [3, 13]]
    for v in vertices:
        assert v['support'] == [i for i, x in enumerate(v['covector']) if Fraction(x)]
        assert all(Fraction(x) <= 0 for x in v['covector'])
    target_vertex = next(i for i, v in enumerate(vertices) if v['support'] == [3, 13])
    saved_faces = {tuple(f['white']): f for f in old['face_ledger']['two_faces']}
    rows = []
    for pair in combinations(sorted(CUSP), 2):
        allowed = set(pair) | {3, 13}
        white = NODES - allowed
        ids = [i for i, v in enumerate(vertices) if set(v['support']) <= allowed]
        sides = [i for i in sorted(allowed)
                 if sum(not Fraction(vertices[j]['covector'][i]) for j in ids) >= 2]
        assert target_vertex in ids and len(ids) == len(sides) in (3, 4)
        frozen = saved_faces[tuple(sorted(white))]
        assert frozen['vertices'] == ids and frozen['side'] == sides
        corners = []
        for j in ids:
            v = vertices[j]
            marks = sorted(set(sides) - set(v['support']))
            assert len(marks) == 2
            actual_star = white | set(marks)
            comp, ds = component_and_distances(actual_star, marks[0])
            distance = ds.get(marks[1])
            if j == target_vertex:
                assert set(marks) == set(pair) and comp == CUSP
                payment = comparison_payment = None
            elif v['ideal']:
                assert v['support'] == [9]
                assert distance is not None and comp in (LEFT, RIGHT)
                assert len(comp) == 9
                payment = comparison_payment = 1  # Fixed whole-factor<=16 unit baseline.
            else:
                comparison_payment = int(distance is not None and distance <= 9)
                payment = int(distance is not None and (distance <= 9 or len(comp) <= 17))
            corners.append(dict(vertex=j, vertex_support=v['support'], marks=marks,
                ideal=v['ideal'], component=sorted(comp), distance=distance,
                fixed_payment=payment, distance9_only_payment=comparison_payment,
                target_D16=j == target_vertex))
        demand = 2 if len(ids) == 3 else 1
        other = sum(c['fixed_payment'] for c in corners if not c['target_D16'])
        deficit = max(0, demand - other)
        assert deficit in (0, 1)
        comparison_other = sum(c['distance9_only_payment'] for c in corners if not c['target_D16'])
        comparison_deficit = max(0, demand - comparison_other)
        assert deficit <= comparison_deficit <= 1
        target = next(c for c in corners if c['target_D16'])
        all_quad_carriers = len(ids) == 4 and LEFT | RIGHT | {9} == white | set(sides)
        if len(ids) == 4:
            assert all_quad_carriers and 0 not in ids
            assert frozen['small_carrier_factors'] == [0, 1]
        rows.append(dict(marked_pair=list(pair), polygon=len(ids), white=sorted(white),
            side=sides, vertices=ids, corners=corners, required_face_payment=demand,
            other_corner_payment=other, D16_payment_deficit=deficit,
            distance9_only_other_payment=comparison_other,
            distance9_only_D16_deficit=comparison_deficit,
            D16_distance=target['distance'], contains_E8_apex=0 in ids,
            E8_factors_in_original_star=[0, 1] if all_quad_carriers else
                frozen['small_carrier_factors']))
    triangles = [r for r in rows if r['polygon'] == 3]
    quads = [r for r in rows if r['polygon'] == 4]
    required = [r for r in rows if r['D16_payment_deficit']]
    rt = [r for r in triangles if r['D16_payment_deficit']]
    rq = [r for r in quads if r['D16_payment_deficit']]
    long = [r for r in rq if r['D16_distance'] > 9]
    cq = [r for r in quads if r['distance9_only_D16_deficit']]
    ct = [r for r in triangles if r['distance9_only_D16_deficit']]
    clong = [r for r in cq if r['D16_distance'] > 9]
    assert (len(rows), len(triangles), len(quads), len(rt), len(rq), len(long)) == (136, 72, 64, 60, 49, 15)
    assert (len(ct), len(cq), len(clong)) == (60, 62, 26)
    assert sum(r['distance9_only_D16_deficit'] != r['D16_payment_deficit'] for r in rows) == 13
    assert histogram(cq, lambda r: r['D16_distance']) == {
        '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7,
        '9': 8, '10': 7, '11': 6, '12': 5, '13': 4, '14': 4}
    # These generators give every automorphism: two independent fork-leaf
    # interchanges and reversal of the unique double-fork spine.
    def identity(): return {i: i for i in CUSP}
    a, b, rev = identity(), identity(), identity()
    a[1], a[2] = 2, 1
    b[11], b[12] = 12, 11
    for i in range(9):
        if i != 3: rev[i], rev[i+10] = i+10, i
    generators = [a, b, rev]
    for g in generators:
        assert set(g.values()) == CUSP
        assert all((j in GRAPH[i]) == (g[j] in GRAPH[g[i]]) for i in CUSP for j in CUSP)
    ids = sorted(CUSP)
    group = {tuple(ids)}
    pending = [identity()]
    while pending:
        p = pending.pop()
        for g in generators:
            new = {i: g[p[i]] for i in CUSP}
            key = tuple(new[i] for i in ids)
            if key not in group: group.add(key); pending.append(new)
    assert len(group) == 8
    actions = [dict(zip(ids, g)) for g in sorted(group)]
    mandatory = {tuple(r['marked_pair']) for r in required}
    closure = {}
    for pair in combinations(ids, 2):
        witnesses = []
        for index, g in enumerate(actions):
            image = tuple(sorted(g[i] for i in pair))
            if image in mandatory: witnesses.append(dict(automorphism=index, mandatory_image=list(image)))
        assert witnesses
        closure[pair] = witnesses[0]
    paths = [Path(__file__)]
    for stem in (OLD, MIXED):
        paths += [ROOT/f'scripts/{stem}.py', ROOT/f'results/{stem}.json', ROOT/f'results/{stem}.md']
    return dict(status='PASS_KNOWN_H17_D16_PAYMENT_BOTTLENECK',
        sources={str(p.relative_to(ROOT)): sha(p) for p in paths},
        hypotheses=dict(finite_corner_payment='unit if marked distance<=9 OR connected marked component order<=17, zero otherwise',
            comparison_finite_corner_payment='unit if marked distance<=9, zero otherwise',
            other_ideal_corner_payment='unit on the actual nine-node E8 varying factor',
            D16_weights='arbitrary nonnegative real weights; orbit lower bound assumes intrinsic diagram invariance'),
        cusp=dict(type='affine_D16', diagram_order=17, nodes=ids,
            other_cusp='affine_E8 x affine_E8', known_finite_volume=True),
        faces=rows, summary=dict(triangles=72, quads=64, mandatory_triangle_pairs=60,
            mandatory_quad_pairs=49, mandatory_long_quad_pairs=15,
            context_specific_mandatory_pairs=109, context_specific_optional_pairs=27,
            mandatory_triangles_with_other_ideal_corner=sum(r['contains_E8_apex'] for r in rt),
            mandatory_quad_distance_histogram=histogram(rq, lambda r:r['D16_distance'])),
        distance9_only_comparison=dict(mandatory_triangle_pairs=60, mandatory_quad_pairs=62,
            mandatory_long_quad_pairs=26, mandatory_pairs=122, optional_pairs=14,
            quads_changed_by_retained_small17_rule=13,
            mandatory_quad_distance_histogram=histogram(cq, lambda r:r['D16_distance']),
            orbit_closure_pairs=136,
            closure_reason='Every primary mandatory pair is also comparison-mandatory, and the saved primary orbit demands already cover all136 pairs.'),
        automorphisms=[{str(i): g[i] for i in ids} for g in actions],
        orbit_demands=[dict(pair=list(p), **w) for p, w in closure.items()],
        invariant_D16_cost_lower_bound=136, normalized_cost='8',
        cross_cusp_scope='All64 quads have both E8 factors in their original star, at a nonincident apex. No bounded cross-cusp multiplicity or global inequality is inferred.',
        scope='A reproducible payment obstruction in an already frozen known finite-volume chamber. No new chamber, dimension bound, or universal lower bound for context-dependent weights.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    d = json.loads(json.dumps(run())); p = ROOT/f'results/{TAG}.json'
    if args.verify: assert d == json.loads(p.read_text())
    else: p.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], '136 faces / 109 mandatory / 122 comparison / 8 automorphisms / 136 orbit demands', flush=True)
